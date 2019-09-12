#!/usr/bin/bash

kubectl exec -it --container=main $1 -- /bin/bash
