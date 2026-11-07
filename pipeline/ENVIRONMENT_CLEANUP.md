# Environment and evidence cleanup work

- The earlier closeout removed 46 native Python environments, 19 native R libraries and a staged Conda installation, with about 69 GiB of observed additional F: free space at that time.
- The science WSL distribution’s 114 scientific Conda environments, recognized R packages and caches were removed. About 204 GiB was freed inside that filesystem; trimming succeeded, while VHD compaction was not performed because it required administrator access.
- Two stopped scientific audit containers and four audit images were removed after unique container outputs had initially been preserved. Unrelated Docker services and volumes were retained.
- The C: Basilisk and MS2PIP caches and standalone Anaconda and Miniconda installations were removed. Their earlier observed free-space increases were approximately 2.83 GiB and a further 8.07 GiB respectively.
- Sam subsequently ended evidence retention. Detailed package inventories, uninstall logs, recovered-output archives, audit results, evaluation trajectories and duplicate Git bundles were deleted in favor of prose accounts.
- Audit data, snapshots, logs and derived results were removed from the retired working areas. Historical setup source and remaining standalone tools are separate from the discarded execution trails.
- Shared runtime software, unrelated Petri Dish work, source repositories and ordinary Git history were not erased by this cleanup.
