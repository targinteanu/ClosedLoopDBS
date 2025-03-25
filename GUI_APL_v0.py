import paramiko
from datetime import datetime
import os
import shutil
import time
from ar_model_fitting import ar_from_csv

# Function to copy only the most recent file from the device to the specified folder
def copy_files_from_device(remote_path, local_path):
    sftp = client.open_sftp()
    try:
        file_list = sftp.listdir_attr(remote_path)
        if not file_list:
            print("UI: No files found in the remote directory.")
            return None
        
        # Find the most recent file based on modification time
        most_recent_file = max(file_list, key=lambda x: x.st_mtime)
        remote_file = os.path.join(remote_path, most_recent_file.filename)
        local_file = os.path.join(local_path, most_recent_file.filename)
        
        sftp.get(remote_file, local_file)
        print(f"UI: Copied {remote_file} to {local_file}")
        return most_recent_file.filename
    finally:
        sftp.close()

# Get the folder location from the user
local_folder = input("Enter the folder location on this computer to copy files to: ")

# Ensure the local folder exists
if not os.path.exists(local_folder):
    os.makedirs(local_folder)

# Specify the remote folder path
remote_folder = '/home/jhuadmin/app/neuromod_software/output/'
remote_yaml = '/home/jhuadmin/app/neuromod_software/configs/realtime_analog_recording.yaml'

def run_neuromod(myShell):
    myShell.send("./run_neuro_modulation.sh realtime_analog_recording.yaml\n")
    time.sleep(2)  # Adjust this delay if needed
    output = myShell.recv(1024).decode()
    print(output)

def stop_neuromod(myShell):
    myShell.send("\x03")  # Send Ctrl+C
    time.sleep(2)  # Allow time for program to terminate
    output = myShell.recv(1024).decode()
    print(output)

# function to modify threshold 
def modify_threshold(myClient, new_thresh, remote_path, local_path):
    print("UI: Modifying threshold value in the YAML file")

    # round new_threshold to the nearest integer
    new_thresh = round(new_thresh)

    # Download the YAML file
    sftp = myClient.open_sftp()
    local_yaml = os.path.join(local_path, os.path.basename(remote_path))
    sftp.get(remote_path, local_yaml)

    # Read the YAML file
    with open(local_yaml, 'r') as file:
        lines = file.readlines()

    # Modify line 52
    lines[51] = f"    thresholdSetting                 : {new_thresh} # FFT Magnitude threshold. Unitless. Must be tuned using threshold results during recording (don't forget the THRESHOLD_DETERMINATION_BLOCK_EXPONENT value [i.e the binary left shift] in the results when trying to use the output to determine this number)\n"
    currentthresh = lines[51]

    # Write the modified content back to the YAML file
    with open(local_yaml, 'w') as file:
        file.writelines(lines)

    # Upload the modified YAML file back to the remote server
    sftp.put(local_yaml, remote_path)
    sftp.close()

    return currentthresh

# function to modify AR coefficients 
def modify_ar_coefficients(myClient, NewArCoeffs, remote_path, local_path):
    print("UI: Modifying AR coefficients in the YAML file")

    # Download the YAML file
    sftp = myClient.open_sftp()
    local_yaml = os.path.join(local_path, os.path.basename(remote_path))
    sftp.get(remote_path, local_yaml)

    # Read the YAML file
    with open(local_yaml, 'r') as file:
        lines = file.readlines()

    # Modify the YAML content
    for i in range(40, 50):
        #print(lines[i])
        lines[i] = f"    arCoeff{i-40} : {NewArCoeffs[i-40]}\n"

    # read the curent threshold value
    currentthresh = lines[51]

    # Write the modified content back to the YAML file
    with open(local_yaml, 'w') as file:
        file.writelines(lines)

    # Upload the modified YAML file back to the remote server
    sftp.put(local_yaml, remote_path)
    sftp.close()

    return currentthresh

# Device credentials
hostname = '192.168.50.70'
username = 'jhuadmin'
password = '123qaz!@#QAZ' 

# Create an SSH client
client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

try:
    # Connect to the device
    print("UI: attempting to connect to device")
    client.connect(hostname, username=username, password=password)
    print("UI: Connected successfully to", hostname)
    time.sleep(1)
    # Start an interactive shell session
    shell = client.invoke_shell()
    time.sleep(5)
    output = shell.recv(1024).decode()
    print(output)

    # Get the current datetime of this computer
    current_datetime = datetime.now().strftime('%Y%m%d %H:%M:%S')
    # Set the date of the device
    print("UI: setting date and time of device")
    shell.send('date \n')
    time.sleep(2)
    output = shell.recv(1024).decode()
    print(output)
    command = f'sudo -S date --set="{current_datetime}" \n'
    shell.send(command)
    time.sleep(2)
    shell.send(password + '\n')
    time.sleep(1)
    output = shell.recv(1024).decode()
    print(output)
    shell.send('date \n')
    time.sleep(5)
    output = shell.recv(1024).decode()
    print(output)
    
    # run the neuro modulation script for initial data
    print("UI: running neuro modulation for initial data")
    shell.send("cd app/neuromod_software\n")
    time.sleep(1)
    run_neuromod(shell)
    time.sleep(10)
    # Stop the neuro modulation script
    print("UI: stopping neuro modulation")
    stop_neuromod(shell)
    time.sleep(5)
    # Copy the most recent file from the device to the specified folder
    print("UI: copying CSV file")
    csvfilename = copy_files_from_device(remote_folder, local_folder)
    time.sleep(5)
    print("UI: getting AR coefficients...")
    arcoeffs = ar_from_csv(local_folder, csvfilename)
    print(arcoeffs)
    curthreshinfo = modify_ar_coefficients(client, arcoeffs, remote_yaml, local_folder)
    # Print line 52 of the YAML file, i.e. current threshold value
    print("Line 52:", curthreshinfo)

    # run the neuro modulation script for additional data
    print("UI: resuming neuro modulation")
    run_neuromod(shell)
    time.sleep(1)

    while True:
        user_input = input("Enter 's' to save data, 'r' to re-set AR coefficients, 't' to adjust threshold, or 'q' to quit: ").strip().lower()
        
        if (user_input == 's') or ((user_input == 'r') or (user_input == 't')):
            # Stop the neuro modulation script
            print("UI: stopping neuro modulation")
            stop_neuromod(shell)
            time.sleep(1)

            if user_input == 'r':
                # modify the AR coefficients 
                print("UI: copying CSV file")
                csvfilename = copy_files_from_device(remote_folder, local_folder)
                print("UI: getting AR coefficients...")
                arcoeffs = ar_from_csv(local_folder, csvfilename)
                print(arcoeffs)
                curthreshinfo = modify_ar_coefficients(client, arcoeffs, remote_yaml, local_folder)

            elif user_input == 't':
                # Print line 52 of the YAML file, i.e. current threshold value
                print("Line 52:", curthreshinfo)
                # Update the threshold value
                new_threshold = input("Enter the new threshold value: ").strip()
                new_threshold = round(float(new_threshold))
                # Modify line 52 of the YAML file
                curthreshinfo = modify_threshold(client, new_threshold, remote_yaml, local_folder)

            # aquire new data with the updated settings 
            print("UI: running neuro modulation for additional data")
            run_neuromod(shell)
            time.sleep(10) 

            # Stop the neuro modulation script
            print("UI: stopping neuro modulation")
            stop_neuromod(shell)
            time.sleep(1)
            # Copy the most recent file for analysis elsewhere 
            print("UI: copying CSV file")
            csvfilename = copy_files_from_device(remote_folder, local_folder)

            # leave the neuro modulation script running 
            print("UI: resuming neuro modulation")
            run_neuromod(shell)
            time.sleep(1)          

        elif user_input == 'q':
            print("UI: Terminating the loop.")
            # Stop the neuro modulation script
            print("UI: stopping neuro modulation")
            stop_neuromod(shell)
            time.sleep(1)
            # Copy the most recent file for analysis elsewhere 
            print("UI: copying CSV file")
            csvfilename = copy_files_from_device(remote_folder, local_folder)
            break
        else:
            print("Invalid input. Please enter 'r', 't', or 'q'.")

finally:
    user_input = input("Power off device? [y/n]: ").strip().lower()
    if user_input == 'y':
        print("UI: powering off device")
        shell.send("sudo poweroff \n")
        time.sleep(2)
        shell.send(password + '\n')
        time.sleep(2)
        output = shell.recv(1024).decode()
        print(output)
    else:
        # Close the connection
        client.close()
        print("UI: Connection closed")