#!/bin/bash

handle_sigterm() {
   echo "SIGTERM received"
   criu dump -t $pid -D /home/baccala/helium/bertini/checkpoint
   exit 0
}

trap handle_sigterm SIGTERM

cd /home/baccala/helium/bertini/run
/BertiniLinux64_v1.6/bertini &
pid=$$
wait
