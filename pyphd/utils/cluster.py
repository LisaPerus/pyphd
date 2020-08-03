# Create and run pbs files on cluster.
# Copyright (C) 2018  Lisa Perus
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# long with this program.  If not, see <https://www.gnu.org/licenses/>.

# System import
import os
import subprocess
import math
import time
import multiprocessing


def create_pbs_cmd(
        cmd,
        pbs_outdir,
        kwargs_type,
        job_name,
        job_number,
        cluster_logdir,
        nb_cpus="1",
        cluster_memory="1",
        cluster_walltime="2",
        cluster_queue=None,
        **kwargs):
    """ Create a pbs command files.

    Parameters:
    -----------
    cmd: list of str
        bash command to be written in pbs file.
    pbs_outdir:
        directory where to store pbs files.
    kwargs_type: dict
        ditc containing for each arg in kwargs its type.
        An argument can be 'short' ('-' will be added to the command)
        or long ('--' will be added to the command).
    kwargs:
      the script parameters: iterative kwargs must contain a list of elements
        and must all have the same length, non-iterative kwargs will be
        replicated.
      If element is None only argname is added. E.g : "T" : [None, None]
      -> -T .
    job_name: str
        job name.
    job_number: int
        job number.
    cluster_logdir: str
        an existing path where the cluster error and output files will be
        stored. This folder must be empty.
    nb_cpus: str (optional, default "1")
        the number of cpus to be used.
    cluster_memory: str (optional, default "1")
        the memory allocated to each job submitted on a cluster (in GB).
    cluster_walltime: str (optional, default "2")
        the walltime used for each job submitted on the cluster (in hours).
    cluster_queue: str (optional, default None)
        the name of the queue where the jobs will be submited.

    Returns
    -------
    pbs_file: str
        pbs file containing command to run script.
    cmd: list of str
        final command.
    log_file: str
        pbs output log file
    err_file: str
        pbs error log file
    """

    if not os.path.isdir(cluster_logdir):
        print("Creating log dir : {0}...".format(cluster_logdir))
        os.mkdir(logdir)

    PBS_TEMPLATE = """#!/bin/bash
#PBS -l mem={memory}gb,nodes=1:ppn={threads},walltime={walltime}:00:00
#PBS -N {name}
#PBS -e {errfile}
#PBS -o {logfile}"""

    # Set cluster parameters
    pbs_script = PBS_TEMPLATE.replace("{memory}", cluster_memory)
    pbs_script = pbs_script.replace("{threads}", nb_cpus)
    pbs_script = pbs_script.replace("{walltime}", cluster_walltime)
    pbs_script = pbs_script.replace("{name}", job_name)
    log_file = os.path.join(cluster_logdir, "{0}_{1}.out.txt".format(
        job_name, job_number))
    err_file = os.path.join(cluster_logdir, "{0}_{1}.err.txt".format(
        job_name, job_number))
    pbs_script = pbs_script.replace("{logfile}", log_file)
    pbs_script = pbs_script.replace("{errfile}", err_file)
    if cluster_queue is not None:
        pbs_script = pbs_script + "\n#PBS -q {0}".format(cluster_queue)

    # Set command
    for arg, arg_value in kwargs.items():
        if kwargs_type[arg] == "short":
            type_arg = "-"
        elif kwargs_type[arg] == "long":
            type_arg = "--"
        else:
            raise ValueError("Unknown argtype : {0}".format(kwargs_type[arg]))
        cmd += [type_arg + arg]
        if arg_value[job_number] is not None:
            cmd += [str(arg_value[job_number])]

    pbs_script = pbs_script + "\n" + " ".join(cmd)

    # Write PBS file
    pbs_file = os.path.join(
        pbs_outdir, "{0}_{1}.pbs".format(job_name, job_number))
    with open(pbs_file, "wt") as open_file:
        open_file.write(pbs_script)

    return pbs_file, cmd, log_file, err_file


def run_qsub(pbs_file):
    """Python wrapper to qsub command.

    Parameters:
    -----------
    pbs_file: str
        path to pbs file to run. All parameters are supposed to have been set
        inside the .pbs file.

    Returns
    -------
    stdout:
        qsub output
    stderr:
        qsub error msg
    error_code: int
        qsub command error code
    """
    cmd = ["qsub", pbs_file]
    cmd = " ".join(cmd)
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    error_code = process.returncode

    return stdout, stderr, error_code


def run_jobs_batch(pbs_files, error_files, cmds, user, queue,
                   nb_jobs_batch=100, logfile=None):
    """Function to run set of jobs

    Divide jobs into batch and run them.
    Jobs are only run is jobs from user are not running yet.
    If user has already jobs running, waits until jobs are finished.

    Parameters:
    -----------
    pbs_files: list of str
        list of pbs files
    error_files: list of str
        list of error files outputs
    cmds: list of str
        list of commands in pbs files
    user: str
        username
    queue: str
        queue name
    nb_jobs_batch: int
        number of jobs per batch
    logfile: str
        path to logfile to write into (optional)
    """

    # Divide jobs into batches
    batches = []
    new_batch = []
    cpt = 0
    for idx, fid in enumerate(pbs_files):
        if cpt == nb_jobs_batch:
            batches.append(new_batch)
            new_batch = []
            cpt = 0
        new_batch.append(fid)
        cpt += 1

    if (len(pbs_files) < nb_jobs_batch or
            ((len(pbs_files) % nb_jobs_batch) != 0)):
        batches.append(new_batch)

    # Run the batches
    cpt = 0
    for batch in batches:
        while True:

            # > Check if user has jobs running
            jobs_running = get_user_queue_jobs(user, queue)
            if len(jobs_running) > 0:
                time.sleep(5)

            # > Else run jobs
            else:
                jobs_ids = []
                qsub_error_msgs = []
                qsub_error_codes = []
                with multiprocessing.Pool() as pool:
                    all_results = pool.map(run_qsub, batch)
                for result in all_results:
                    job_id = result[0].decode("utf-8").strip("\n")
                    qsub_error_msg = result[1].decode("utf-8")
                    qsub_error_code = result[2]
                    jobs_ids.append(job_id)
                    qsub_error_msgs.append(qsub_error_msg)
                    qsub_error_codes.append(qsub_error_code)

                # >> Wait until jobs are finished to run logs
                while True:
                    if check_jobs_finished_running(jobs_ids):

                        # If log file, write logs
                        if logfile is not None:
                            with open(logfile, "at") as open_file:
                                for idx, pbs_file in enumerate(batch):
                                    line = "cluster_" + os.path.basename(
                                        pbs_file)
                                    line += "_cmd = "
                                    line += str(cmds[cpt])
                                    line += "]\n"
                                    line += "cluster_" + os.path.basename(
                                        pbs_file)
                                    if qsub_error_codes[idx] != 0:
                                        error_code = qsub_error_codes[idx]
                                        error_msg = qsub_error_msgs[idx]
                                    else:
                                        qsub_err_file = error_files[cpt]
                                        with open(
                                                qsub_err_file, "rt") as of:
                                            error_lines = of.readlines()
                                        if len(error_lines) == 1:
                                            if len(error_lines[0].strip(
                                                    " ").strip("\n")) == 0:
                                                error_code = 0
                                                error_msg = ""
                                            else:
                                                error_code = 1
                                                error_msg = " - ".join(
                                                    error_lines)
                                        else:
                                            error_code = 0
                                            error_msg = ""
                                    line += "_exitcode = " + str(error_code)
                                    line += "\n"
                                    line += "cluster_" + os.path.basename(
                                        pbs_file)
                                    line += "_error = " + "[" + error_msg + "]"
                                    open_file.write(line)
                                    open_file.write("\n")
                                    print(line)
                                    cpt += 1
                        break
                break


def get_user_queue_jobs(user, queue):
    """Function to get jobs IDs from a queue.

    Parameters:
    -----------
    user: str
        user name on cluster.
    queue: str
        queue name.

    Returns
    -------
    jobs_ids: list of str
    """
    cmd = "qstat | grep {0} | grep {1}".format(user, queue)
    cmd += " | awk '{print $1}'"
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    jobs_ids = stdout.decode("utf-8").split("\n")
    jobs_ids = [x for x in jobs_ids if len(x) != 0]
    return jobs_ids


def check_jobs_finished_running(jobs):
    """Function to check if list of jobs has finished running.

    Parameters:
    -----------
    jobs: list of str
        list of jobs.

    Returns
    -------
    jobs_run: bool
        True if all jobs have been run.
    """
    all_job_run = []
    for job in jobs:
        cmd = "qstat -x {0}".format(job)
        cmd += " | awk '{print $5}'"
        process = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        job_status = stdout.decode("utf-8").split("\n")[-1].strip(" ")
        if job_status == "F":
            all_job_run.append(True)
        else:
            all_job_run.append(False)
    jobs_run = True
    for status in all_job_run:
        if not status:
            jobs_run = False
    return jobs_run
