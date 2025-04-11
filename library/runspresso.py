from subprocess import getoutput
import numpy as np
import subprocess
import time
import re

class runspresso:

  def __init__(self):
    self.a = 0

  def create_script(self, job_id='nan', nNodes=1):
    script = f"""#!/bin/bash
mpirun -np NNN pw.x -in {job_id}.in > {job_id}.out
      """
    script = script.replace('NNN', f'{nNodes}')
    return script
    
  # Function to run a command in a detached screen session
  def run_job_in_screen(self, job_id='nan', nNodes=1):
    script = self.create_script(job_id, nNodes)
    script_file = f"run_{job_id}.sh"
    with open(script_file, "w") as file:
      file.write(script)

    subprocess.run(f"chmod +x {script_file}", shell=True)
    subprocess.run(f"screen -dmS job_{job_id} nice -20 ./{script_file}", shell=True)
    pid = getoutput(f"screen -ls | grep 'job_{job_id}' | awk '{{print $1}}' | cut -d. -f1")
    print("Running on screen:", pid)
    return(pid)

  def check_screen_sessions(self, job_id='nan'):
    result = subprocess.run("screen -ls", shell=True, capture_output=True, text=True)
    active_screens = result.stdout
    if f"job_{job_id}" in active_screens:
      return True
    else:
        return False

  def check_if_all_finished(self, n_indv=1):
    while True:
      finished = [True for _ in range(n_indv)]
      for i in range(n_indv):
        finished[i] = self.check_screen_sessions(i+1)
      if not any(finished):
        print("All jobs are finished.")
        print(finished)
        break  # Exit the loop if all jobs are finished
      else:
        print("Yucas. Jobs active:", finished)
      time.sleep(10)

  def fetch_total_energy(self, infile):
    total_energy = None
    try:
      with open(infile, 'r') as file:
        for line in file:
          if '!' in line and 'total energy' in line:            
            match = re.search(r'total energy\s+=\s+([-0-9.]+)\s+Ry', line)
            if match:
              total_energy = float(match.group(1))
              break
    except FileNotFoundError:
      print(f"File Not Found:", infile)
    except Exception as e:
      print(f"Error While Reading File {file}: {e}")
    return total_energy

  def get_total_energies(self, indv_key='nan'):
    ener = self.fetch_total_energy(f'{indv_key}.out')
    if ener is not None:
      return ener
    else:
      return -1