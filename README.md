# Purpose:

This is a template for a website that runs Python on the backend

# Install Modules:

You will need to install the following Python modules:

    pip install Flask

# Terminate Server:

Leaving the process running can be a significant drain on your battery, so its best to turn things off after you are done.

- MacOS:
    List open files on port 5000:

        lsof -i :5000

    Use the process ID to terminate the Python:

        kill 

- Windows

    Inspect network connections on port 5000

        netstat -ano | findstr :5000
    
    - Command output:

            TCP    127.0.0.1:5000         0.0.0.0:0              LISTENING       20352

    Lists all running tasks with matching process ID:

        tasklist /FI "PID eq 20352"

    - Command output:
  
            Image Name                     PID Session Name        Session#    Mem Usage
            ========================= ======== ================ =========== ============
            python.exe                   20352 Console                    3     39,732 K

    Kill the process:

        taskkill /PID 20352 /F

    - Command output:

            SUCCESS: The process with PID 20352 has been terminated.
