# system tests for toast CLI functions

import subprocess
import os
import filecmp
import pytest

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
    
def test_design_400_25nn_userinput(setup_environment):
    # Resolve path to test_df.csv relative to this test file
    test_dir = os.path.dirname(__file__)
    user_data = os.path.join(test_dir, "test_df.csv")

    command = (
        f"toast design "
        f"-op {setup_environment} "
        f"-a 400 "
        f"-sn 1 "
        f"-sg rpoB,katG "
        f"-nn 2 "
        f"-p 150 "
        f"-ud {user_data}"
    )
    returncode, stdout, stderr = run_command(command)

    print("STDOUT:", stdout.decode())
    print("STDERR:", stderr.decode())
    assert returncode == 0

    # Add more assertions to check the expected output
def test_design_400_25nn(setup_environment):
    command = (
        f"toast design "
        f"-op {setup_environment} "
        f"-a 400 "
        f"-sn 1 "
        f"-sg rpoB,katG "
        f"-nn 2 "
        f"-p 150"

    )
    returncode, stdout, stderr = run_command(command)

    assert returncode == 0, stderr.decode()


def test_design_400_25sp(setup_environment):
    # Example for one variant of the design function
    command = f"toast design -op {setup_environment} -a 400 -sn 2 -sg rpoB,katG -nn 1"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output
    
def test_design_1000_25sp(setup_environment):
    # Example for one variant of the design function
    command = f"toast design -op {setup_environment} -a 1000 -sn 1 -sg rpoB,katG -nn 2"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output

def test_design_1000_25sp_sc(setup_environment):
    # Example for one variant of the design function
    command = f"toast design -op {setup_environment} -a 1000 -sn 2 -sg rpoB,katG -nn 1 -sc"
    returncode, stdout, stderr = run_command(command)
    assert returncode == 0
    # Add more assertions to check the expected output
    

def test_dummy():
    assert True


#non design function tests

# def test_plotting_function(setup_environment):
#     # Prepare input files as required for the plotting function
#     amplicon_file = setup_environment / "Primer_design-accepted_primers-26-400.csv"
#     reference_file = "../db/reference_design.csv"
#     command = f"toast plotting -ap {amplicon_file} -rp {reference_file} -op {setup_environment} -r 400"
#     returncode, stdout, stderr = run_command(command)
#     assert returncode == 0
# from pathlib import Path

# def test_plotting_function(setup_environment):
#     # Step 1: generate required design output
#     design_cmd = (
#         f"toast design "
#         f"-op {setup_environment} "
#         f"-a 400 "
#         f"-sn 1 "
#         f"-sg rpoB,katG "
#         f"-nn 2 "
#         f"-p 150"
#     )
#     rc, stdout, stderr = run_command(design_cmd)
#     assert rc == 0, stderr.decode()

#     # Step 2: plotting
#     amplicon_file = setup_environment / "Primer_design-accepted_primers-3-400.csv"
#     assert amplicon_file.exists(), f"Missing amplicon file: {amplicon_file}"

#     project_root = Path(__file__).resolve().parents[1]
#     reference_file = project_root / "toast_amplicon" / "db" / "reference_design.csv"
#     assert reference_file.exists(), f"Missing reference file: {reference_file}"

#     plot_cmd = (
#         f"toast plotting "
#         f"-ap {amplicon_file} "
#         f"-rp {reference_file} "
#         f"-op {setup_environment} "
#         f"-r 400"
#     )

#     rc, stdout, stderr = run_command(plot_cmd)
#     print("STDOUT:", stdout.decode())
#     print("STDERR:", stderr.decode())

#     assert rc == 0

    # Add more assertions to check the expected output
    
# def test_amplicon_no(setup_environment):
#     command = f"toast amplicon_no -a 400 -op {setup_environment} -g"
#     returncode, stdout, stderr = run_command(command)
#     assert returncode == 0
#     # Add more assertions to check the expected output