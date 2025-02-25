import paramiko
from datetime import datetime
import time

# Device credentials
hostname = '192.168.50.70'
username = 'jhuadmin'
password = '123qaz!@#QAZ' 

# Establish SSH connection
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

try:
    ssh.connect(hostname, username=username, password=password)
    
    # Start an interactive shell session
    shell = ssh.invoke_shell()

    # Wait for shell to be ready
    time.sleep(1)

    # Get the current datetime of this computer
    current_datetime = datetime.now().strftime('%Y%m%d %H:%M:%S')

    # Set the date of the device
    command = f'sudo -S date --set="{current_datetime}"'
    shell.send(command + '\n')
    time.sleep(1)
    shell.send(password + '\n')
    time.sleep(1)
    output = shell.recv(1024).decode()
    print(output)
    
    # Send commands
    shell.send("cd app/neuromod_software\n")
    time.sleep(1)
    shell.send("./run_neuro_modulation.sh realtime_analog_recording.yaml\n")

    # Read output
    time.sleep(2)  # Adjust this delay if needed
    output = shell.recv(1024).decode()
    print(output)

    print("Program is running. Press Enter to stop...")
    
    # Keep reading output
    while True:
        if shell.recv_ready():
            output = shell.recv(1024).decode()
            print(output, end="")

        # Stop when the user presses Enter
        if input().strip() == "":
            print("\nStopping the program...")
            shell.send("\x03")  # Send Ctrl+C
            break

    time.sleep(2)  # Allow time for program to terminate
    output = shell.recv(1024).decode()
    print(output)

    ssh.close()

except KeyboardInterrupt:
    print("\nInterrupted.")
    ssh.close()