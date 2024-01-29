import pytest
import subprocess
import os
import filecmp

@pytest.fixture(scope='module')
def setup_environment(tmp_path_factory):
    # Setup for test environment
    tmp_dir = tmp_path_factory.mktemp('data')
    output_path = tmp_dir / 'Amplicon_design_output'
    os.makedirs(output_path, exist_ok=True)
    return output_path

def run_command(command):
    # Helper function to run CLI commands
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return process.returncode, stdout, stderr

def test_amplicon_no(setup_environment):
    command = f"toast amplicon_no -a 400 -op {setup_environment} -g"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output
    
def test_design_400_25nn_userinput(setup_environment):
    # Example for one variant of the design function
    command = f"toast design -op {setup_environment} -a 400 -sn 1 -sg rpoB,katG -nn 25 -ud /mnt/storage10/lwang/Projects/TOAST/cache/test_df.csv"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output
def test_design_400_25nn(setup_environment):
    # Example for one variant of the design function
    command = f"toast design -op {setup_environment} -a 400 -sn 1 -sg rpoB,katG -nn 25"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output

def test_design_400_25nn(setup_environment):
    # Example for one variant of the design function
    command = f"toast design -op {setup_environment} -a 400 -sn 1 -sg rpoB,katG -nn 25"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output

def test_design_400_25sp(setup_environment):
    # Example for one variant of the design function
    command = f"toast design -op {setup_environment} -a 400 -sn 25 -sg rpoB,katG -nn 1"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output
    
def test_design_1000_25sp(setup_environment):
    # Example for one variant of the design function
    command = f"toast design -op {setup_environment} -a 1000 -sn 1 -sg rpoB,katG -nn 25"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output

def test_design_1000_25sp_sc(setup_environment):
    # Example for one variant of the design function
    command = f"toast design -op {setup_environment} -a 1000 -sn 25 -sg rpoB,katG -nn 1 -sc"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output
    

def test_plotting_function(setup_environment):
    # Prepare input files as required for the plotting function
    amplicon_file = setup_environment / "Primer_design-accepted_primers-26-400.csv"
    reference_file = "../db/reference_design.csv"
    command = f"toast plotting -ap {amplicon_file} -rp {reference_file} -op {setup_environment} -r 400"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output

# Additional test functions for other variants of the commands...
