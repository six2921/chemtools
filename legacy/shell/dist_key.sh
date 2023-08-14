#!/bin/bash

cat list_machines | while read mcn;do
    scp authorized_keys $user@$mcn:/home/$user/.ssh/
done

