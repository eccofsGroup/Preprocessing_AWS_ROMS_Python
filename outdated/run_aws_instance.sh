#!/bin/bash

aws ec2 start-instances --instance-ids i-0becb62354fcacf73

ssh -i "oceanmodeling.pem" ubuntu@ec2-44-209-181-39.compute-1.amazonaws.com 'python3 mainCodeRoms.py'
scp -i "oceanmodeling.pem" ubuntu@ec2-44-209-181-39.compute-1.amazonaws.com:/home/ubuntu/urlDownloads.txt .
bash urlDownloads.txt

aws ec2 stop-instances --instance-ids i-0becb62354fcacf73
