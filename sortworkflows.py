#!/usr/bin/env python
import argparse
import re

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--workflows', help='The list of workflows')
    args = parser.parse_args()
    return args

def parse_workflows(workflow_file):
    workflows = []
    with open(workflow_file, 'r') as inf:
        for line in inf:
            data = line.split()
            time = data[1]
            days = re.findall(r"(\d+)d", time)
            if len(days) > 0:
                days = int(days[0])
            else:
                days = 0
            hours = re.findall(r"(\d+)h", time)
            if len(hours) > 0:
                hours = int(hours[0])
            else:
                hours = 0
            minutes = re.findall(r"(\d+)m", time)
            if len(minutes) > 0:
                minutes = int(minutes[0])
            else:
                minutes = 0
            total_time = days*1440+hours*60+minutes
            workflows.append((total_time, data))
    return workflows

def main():
    args = parseargs()
    workflows = parse_workflows(args.workflows)
    workflows = sorted(workflows)
    for wf in workflows:
        print(wf[1][0], wf[1][1])

if __name__ == '__main__':
    main()
