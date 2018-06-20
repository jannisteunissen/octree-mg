        subroutine  MPI_INIT(ierr)
        return
        end

        subroutine  MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
        iproc=0
        ierr=0
        return
        end

        subroutine  MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        nproc=1
        ierr=0
        return
        end

        subroutine  MPI_FINALIZE(ierr)
        return
        end

        subroutine  MPI_ALLREDUCE(wrkallred,allred,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        stop 'ALLREDUCE'
        return
        end

        subroutine  MPI_BCAST(psi,nvctr,MPI_DOUBLE_PRECISION,jproc,MPI_COMM_WORLD,ierr)
        stop 'BCAST'
        return
        end

        subroutine  MPI_BARRIER(MPI_COMM_WORLD,ierr)
        return
        end

        subroutine MPI_REDUCE()
        return
        end

        subroutine  MPI_ALLGatherV()
        stop 'ALLGATHERV'
        return
        end

        subroutine  MPI_ALLGather()
        stop 'ALLGATHER'
        return
        end

        subroutine  MPI_ALLTOALL()
        stop 'ALLTOALL'
        return
        end

        subroutine  MPI_REDUCE_SCATTER()
        stop 'REDUCE_SCATTER'
        return
        end

