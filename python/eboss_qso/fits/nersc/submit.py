from __future__ import print_function
import argparse
import os, sys
from jinja2 import Environment, FileSystemLoader, Template
import tempfile
import subprocess
import numpy

# the directory this file lives in
toplevel = os.path.split(os.path.abspath(__file__))[0]

def minutes_to_job_time(minutes):
    h, m = divmod(minutes, 60)
    return "%02d:%02d:00" % (h, m)

def get_nodes_from_cores(cores, host):
    if host == 'cori':
        nodes, extra = divmod(cores, 32)
    elif host == 'edison':
        nodes, extra = divmod(cores, 24)
    else:
        raise ValueError("bad host name '%s'" %host)
    # account for remainder cores
    if extra > 0: nodes += 1
    return nodes

class NERSCManager(object):
    """
    A manager to submit jobs on NERSC.
    """
    @classmethod
    def parse_args(cls):

        parser = cls.get_parser()
        ns, unknown = parser.parse_known_args()

        if ns.subparser_name == 'run':
            return ns

        command = [os.path.basename(sys.executable), sys.argv[0]] + ['run'] + unknown

        # determine the NERSC host
        host = os.environ.get('NERSC_HOST', 'cori')
        if host is None:
             raise RuntimeError("jobs should be executed on NERSC")

        # load the template job script
        jinja_env = Environment(loader=FileSystemLoader(toplevel))
        tpl = jinja_env.get_template('job.template.sh')

        # the configuration to pass to the template
        config = {}
        config['output_file'] = "output/{host}-%j.out.{cores}".format(host=host, cores=ns.cores)
        config['command'] = " ".join(command)
        config['partition'] = ns.partition
        config['time'] = minutes_to_job_time(ns.time)
        config['cores'] = ns.cores
        config['nodes'] = get_nodes_from_cores(ns.cores, host)
        config['haswell_config'] = "#SBATCH -C haswell" if host == 'cori' else ""
        config['job'] = 'eboss-qso-fit'

        # render the template
        rendered = tpl.render(**config)

        # write to temp file and call
        with tempfile.NamedTemporaryFile(mode='w') as ff:

            # echo the job scripts
            print(rendered)

            # write to temp file (and rewind)
            ff.write(rendered)
            ff.seek(0)

            # and call
            subprocess.call(["sbatch", ff.name])

        sys.exit(0)

    @classmethod
    def add_argument(cls, *args, **kwargs):

        parser = cls.get_parser()
        cls._run.add_argument(*args, **kwargs)

    @classmethod
    def get_parser(cls):

        if not hasattr(cls, 'parser'):

            from argparse import ArgumentParser

            parser = ArgumentParser(description='submit a job to NERSC')

            subparsers = parser.add_subparsers(dest='subparser_name')

            h = 'submit a job to NERSC'
            submit = subparsers.add_parser('submit', help=h)

            h = 'the number of nodes to request'
            submit.add_argument('-n', '--cores', type=int, help=h, required=True)

            h = 'the NERSC partition to submit to'
            choices=['debug', 'regular']
            submit.add_argument('-p', '--partition', type=str, choices=choices, default='debug', help=h)

            h = 'the requested amount of time (in minutes)'
            submit.add_argument('-t', '--time', type=int, default=30, help=h)

            cls._run = subparsers.add_parser('run', help='run the job')
            cls.parser = parser

        return cls.parser
