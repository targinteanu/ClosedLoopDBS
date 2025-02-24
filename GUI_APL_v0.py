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
            print("No files found in the remote directory.")
            return None
        
        # Find the most recent file based on modification time
        most_recent_file = max(file_list, key=lambda x: x.st_mtime)
        remote_file = os.path.join(remote_path, most_recent_file.filename)
        local_file = os.path.join(local_path, most_recent_file.filename)
        
        sftp.get(remote_file, local_file)
        print(f"Copied {remote_file} to {local_file}")
        return most_recent_file.filename
    finally:
        sftp.close()

# Get the folder location from the user
local_folder = input("Enter the folder location on this computer to copy files to: ")

# Ensure the local folder exists
if not os.path.exists(local_folder):
    os.makedirs(local_folder)

# Specify the remote folder path
remote_folder = 'jhuadmin@192.168.50.70:/home/jhuadmin/app/neuromod_software/output/'

# function to modify AR coefficients 
def modify_ar_coefficients(myClient, NewArCoeffs):

    # Read the YAML file
    command = 'cat configs/realtime_analog_recording.yaml'
    stdin, stdout, stderr = myClient.exec_command(command)
    yaml_content = stdout.read().decode()

    # Modify the YAML content
    lines = yaml_content.split('\n')
    for i in range(41, 51):
        lines[i] = f"arCoeff{i-41} : {NewArCoeffs[i-41]}"

    new_yaml_content = '\n'.join(lines)

    # Write the modified content back to the YAML file
    command = f'echo "{new_yaml_content}" > configs/realtime_analog_recording.yaml'
    stdin, stdout, stderr = myClient.exec_command(command)
    print(stdout.read().decode())
    print(stderr.read().decode())

# Device credentials
hostname = '192.168.50.70'
username = 'jhuadmin'
password = '123qaz!@#QAZ' 

# Create an SSH client
client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

try:
    # Connect to the device
    client.connect(hostname, username=username, password=password)
    print("Connected successfully to", hostname)
    time.sleep(1)

    # Get the current datetime of this computer
    current_datetime = datetime.now().strftime('%Y%m%d %H:%M:%S')

    # Set the date of the device
    command = f'sudo -S date --set="{current_datetime}"'
    stdin, stdout, stderr = client.exec_command(command)
    stdin.write(password + '\n')
    stdin.flush()
    print(stdout.read().decode())
    print(stderr.read().decode())
    time.sleep(1)
    
    # Change directory to app/neuromod_software
    command = 'cd app/neuromod_software'
    stdin, stdout, stderr = client.exec_command(command)
    print(stdout.read().decode())
    print(stderr.read().decode())

    # run the neuro modulation script for initial data
    command = './run_neuro_modulation.sh realtime_analog_recording.yaml'
    stdin, stdout, stderr = client.exec_command(command)
    print(stdout.read().decode())
    print(stderr.read().decode())
    time.sleep(10)
    # Stop the neuro modulation script
    command = 'pkill -f run_neuro_modulation.sh'
    stdin, stdout, stderr = client.exec_command(command)
    print(stdout.read().decode())
    print(stderr.read().decode())
    time.sleep(1)
    # Copy the most recent file from the device to the specified folder
    csvfilename = copy_files_from_device(remote_folder, local_folder)
    arcoeffs = ar_from_csv(local_folder, csvfilename)
    print(arcoeffs)
    modify_ar_coefficients(client, arcoeffs)
    # Print line 52 of the YAML file, i.e. current threshold value
    command = 'sed -n "52p" configs/realtime_analog_recording.yaml'
    stdin, stdout, stderr = client.exec_command(command)
    print("Line 52:", stdout.read().decode().strip())
    print(stderr.read().decode())

    # run the neuro modulation script for additional data
    command = './run_neuro_modulation.sh realtime_analog_recording.yaml'
    stdin, stdout, stderr = client.exec_command(command)
    print(stdout.read().decode())
    print(stderr.read().decode())
    time.sleep(1)

    while True:
        user_input = input("Enter 'r' to re-set AR coefficients, 't' to adjust threshold, or 'q' to terminate: ").strip().lower()
        
        if (user_input == 'r') or (user_input == 't'):
            # Stop the neuro modulation script
            command = 'pkill -f run_neuro_modulation.sh'
            stdin, stdout, stderr = client.exec_command(command)
            print(stdout.read().decode())
            print(stderr.read().decode())
            time.sleep(1)

            if user_input == 'r':
                # modify the AR coefficients 
                csvfilename = copy_files_from_device(remote_folder, local_folder)
                arcoeffs = ar_from_csv(local_folder, csvfilename)
                print(arcoeffs)
                modify_ar_coefficients(client, arcoeffs)

            elif user_input == 't':
                # Print line 52 of the YAML file, i.e. current threshold value
                command = 'sed -n "52p" configs/realtime_analog_recording.yaml'
                stdin, stdout, stderr = client.exec_command(command)
                print("Line 52:", stdout.read().decode().strip())
                print(stderr.read().decode())
                # Update the threshold value
                new_threshold = input("Enter the new threshold value: ").strip()
                # Modify line 52 of the YAML file
                command = f'sed -i "52s/.*/thresholdSetting                 : {new_threshold}/" configs/realtime_analog_recording.yaml'
                stdin, stdout, stderr = client.exec_command(command)
                print(stdout.read().decode())
                print(stderr.read().decode())

            # aquire new data with the updated settings 
            command = './run_neuro_modulation.sh realtime_analog_recording.yaml'
            stdin, stdout, stderr = client.exec_command(command)
            print(stdout.read().decode())
            print(stderr.read().decode())
            time.sleep(20) 
            # Stop the neuro modulation script
            command = 'pkill -f run_neuro_modulation.sh'
            stdin, stdout, stderr = client.exec_command(command)
            print(stdout.read().decode())
            print(stderr.read().decode())
            time.sleep(1)
            # Copy the most recent file for analysis elsewhere 
            csvfilename = copy_files_from_device(remote_folder, local_folder)

            # leave the neuro modulation script running 
            command = './run_neuro_modulation.sh realtime_analog_recording.yaml'
            stdin, stdout, stderr = client.exec_command(command)
            print(stdout.read().decode())
            print(stderr.read().decode())
            time.sleep(1)          

        elif user_input == 'q':
            print("Terminating the loop.")
            # Stop the neuro modulation script
            command = 'pkill -f run_neuro_modulation.sh'
            stdin, stdout, stderr = client.exec_command(command)
            print(stdout.read().decode())
            print(stderr.read().decode())
            time.sleep(1)
            # Copy the most recent file for analysis elsewhere 
            csvfilename = copy_files_from_device(remote_folder, local_folder)
            break
        else:
            print("Invalid input. Please enter 'r', 't', or 'q'.")

finally:
    # Close the connection
    client.close()
    print("Connection closed")