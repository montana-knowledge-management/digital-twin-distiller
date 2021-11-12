#!/bin/bash
imageName=montanaltd/distiller-wk-2021-1
containerName=distiller-wk-2021-1
echo Delete old container...
docker rm -f $containerName
docker build --no-cache -t $imageName -f Dockerfile  .
echo Run new container...
docker run -d -p 5000:5000 --name $containerName $imageName
