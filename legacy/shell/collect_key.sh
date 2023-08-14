#!/bin/bash
user='siu'
cat list_machines | while read mcn; do
    ssh $user@$mcn cat /home/$user/.ssh/id_rsa.pub
    echo $mcn
done
