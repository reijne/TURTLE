      program expand
      character*256 node,vari
      call getenv("SLURM_JOB_NUM_NODES",vari)
      read(vari,*)numnodes
      call getenv("SLURM_TASKS_PER_NODE",vari)
      ihaak=index(vari,"(")-1
      read(vari(1:ihaak),*)numtask
      do i=1,numnodes
         read(5,*)node
         ispace=index(node," ")-1
         do j=1,numtask
            write(6,'(A)')node(1:ispace)
         enddo
      enddo
      end
      
