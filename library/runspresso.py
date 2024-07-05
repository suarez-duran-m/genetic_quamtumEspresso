import numpy as np
import subprocess
import time
import re

class runspresso:

  def __init__(self):
    self.a = 0

  def create_script(self, job_id='nan', ischild='False'):
    if ischild:
      script = f"""#!/bin/bash
      pw.x -in child{job_id}.in > child_{job_id}.out
      """
    else:
      script = f"""#!/bin/bash
      pw.x -in ind{job_id}.in > {job_id}.out
      """
    return script
    
  # Function to run a command in a detached screen session
  def run_job_in_screen(self, job_id=1, ischild='False'):   
    script = self.create_script(job_id, ischild)
    if ischild:
      script_file = f"run_child_{job_id}.sh"
    else:
      script_file = f"run_ind_{job_id}.sh"
    with open(script_file, "w") as file:
      file.write(script)
        
    subprocess.run(f"chmod +x {script_file}", shell=True)
    subprocess.run(f"screen -dmS job_{job_id} ./{script_file}", shell=True)

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
      print(f"File Not Found.")
    except Exception as e:
      print(f"Error While Reading File {file}: {e}")
    return total_energy

  def get_total_energies(self, n_indv=1, ischild=False):
    total_energies = []
    for ind_id in range(1, n_indv+1):
      if ischild:
        ener = self.fetch_total_energy(f'child_{ind_id}.out')
      else:
        ener = self.fetch_total_energy(f'{ind_id}.out')
      if ener is not None:
        total_energies.append(ener)

    return np.array(total_energies)
