mpiexec	--mca btl_tcp_if_include 10.0.0.0/24 \
        --mca oob_tcp_if_include 10.0.0.0/24 \
	--bind-to none -np 6 --host lemuria:2,node1,node2,node3,node4 $@
