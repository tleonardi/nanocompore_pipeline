all: nanocompore.img

guppy.img: Singularity.guppy
	sudo singularity build guppy.img Singularity.guppy
nanocompore.img: Singularity.nanocompore guppy.img
	sudo singularity build nanocompore.img Singularity.nanocompore
