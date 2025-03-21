import paramiko
from datetime import datetime
import os
import shutil
import time

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
    shell.send("./run_neuro_modulation.sh realtime_analog_recording.yaml\n")
    time.sleep(1)  # Adjust this delay if needed
    output = shell.recv(1024).decode()
    print(output)
    time.sleep(10)

    # Stop the neuro modulation script
    print("UI: stopping neuro modulation")
    shell.send("\x03")  # Send Ctrl+C
    time.sleep(1)  # Allow time for program to terminate
    output = shell.recv(1024).decode()
    print(output)

finally:
    # Close the connection
    client.close()
    print("Connection closed")