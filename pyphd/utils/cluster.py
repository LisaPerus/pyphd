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
    cmd:
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
            cmd += """-{0} {1}
                """.format(arg, arg_value[job_number])
        elif kwargs_type[arg] == "long":
            cmd += """--{0} {1}
                """.format(arg, arg_value[job_number])
    pbs_script = pbs_script + "\n{0}".format(cmd)

    # Write PBS file
    pbs_file = os.path.join(
        pbs_outdir, "{0}_{1}.pbs".format(job_name, job_number))
    with open(pbs_file, "wt") as open_file:
        open_file.write(pbs_script)

    return pbs_file


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
        qsub error code
    """
    cmd = ["qsub", pbs_file]
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    return stdout, stderr
