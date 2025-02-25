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

    # Get the current datetime of this computer
    current_datetime = datetime.now().strftime('%Y%m%d %H:%M:%S')

    # Set the date of the device
    command = f'sudo -S date --set="{current_datetime}"'
    stdin, stdout, stderr = ssh.exec_command(command)
    stdin.write(password + '\n')
    stdin.flush()
    print(stdout.read().decode())
    print(stderr.read().decode())
    time.sleep(1)
    
    # Execute command
    command = "bash -i -c 'cd app/neuromod_software && ./run_neuro_modulation.sh realtime_analog_recording.yaml'"
    stdin, stdout, stderr = ssh.exec_command(command)
    time.sleep(5)

    # Print output in real time
    for line in iter(stdout.readline, ""):
        print(line, end="")

    ssh.close()

except KeyboardInterrupt:
    print("\nInterrupted.")
    ssh.close()