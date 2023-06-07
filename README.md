# SpatialViewR
Prepares data for SpatialView from R. Currently SpatialViewR supports Seurat and SpatialExperiment objects.

### Installlation
SpatialViewR can be installed using 

```r
remotes::install_github("kendziorski-lab/SpatialviewR")
```

### Usage

Please follow
- [A step by step guide to export data from Seurat object](https://kendziorski-lab.github.io/projects/spatialview/SpatialView_Tutorial_Using_Seurat.html)
- [A step by step guide to export data from SpatialExperiment](https://kendziorski-lab.github.io/projects/spatialview/SpatialView_Tutorial_Using_SpatialExperiment.html)


### FAQ

1. Can I change the web port?

Yes. Port number 8878 is default port in the *prepare10x_from_seurat* and *prepare10x_from_SpatialExperiment* functions in [SpatialViewR](https://github.com/kendziorski-lab/SpatialViewR) package. You can pass a different port number in these functions.

2. Once a port is used, I get 'port in use' error on my next run.

Currently, SpatialViewR doesn't track the web process. To avoid the using the previously used port, you can provide a different port number in the *prepare10x_from_seurat* and *prepare10x_from_SpatialExperiment* functions. Alternative, you can kill the current process that uses the port. the following code may be helpful in Linux or MacOS
```
$ lsof -i :8878

#COMMAND   PID      USER   FD   TYPE     DEVICE SIZE/OFF NODE NAME
#Python  55461      XXX    XX  XX.       XX      XX  TCP *:8878 (LISTEN)

$kill 55461
```

