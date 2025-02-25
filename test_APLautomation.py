import paramiko

# Device credentials
hostname = '192.168.50.70'
username = 'jhuadmin'
password = '123qaz!@#QAZ' 

# Command to execute remotely
command = "cd app/neuromod_software && ./run_neuro_modulation.sh realtime_analog_recording.yaml"

# Establish SSH connection
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

try:
    ssh.connect(hostname, username=username, password=password)
    
    # Execute command
    stdin, stdout, stderr = ssh.exec_command(command)

    # Print output in real time
    for line in iter(stdout.readline, ""):
        print(line, end="")

    ssh.close()

except KeyboardInterrupt:
    print("\nInterrupted.")
    ssh.close()