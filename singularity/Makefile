LAST_GUPPY_COMMIT=$(shell git log -n 1 --pretty=format:%h -- Singularity.guppy)
LAST_NANOCOMPORE_COMMIT=$(shell git log -n 1 --pretty=format:%h -- Singularity.nanocompore)
LAST_GUPPY_IMAGE=guppy_$(LAST_GUPPY_COMMIT).img
LAST_NANOCOMPORE_IMAGE=nanocompore_$(LAST_NANOCOMPORE_COMMIT).img

all: nanocompore

.PHONY: nanocompore
nanocompore: $(LAST_NANOCOMPORE_IMAGE)
	
.PHONY: guppy
guppy: $(LAST_GUPPY_IMAGE)

$(LAST_GUPPY_IMAGE): Singularity.guppy
	git diff-index --quiet HEAD -- Singularity.guppy || ( echo "Commit the changes to Singularity.guppy"; exit 1)
	sudo singularity build $(LAST_GUPPY_IMAGE) Singularity.guppy
	ln -s $(LAST_GUPPY_IMAGE) guppy_latest.img

$(LAST_NANOCOMPORE_IMAGE): Singularity.nanocompore $(LAST_GUPPY_IMAGE)
	git diff-index --quiet HEAD -- Singularity.nanocompore || ( echo "Commit the changes to Singularity.nanocompore"; exit 1)
	sudo singularity build --writable $(LAST_NANOCOMPORE_IMAGE) Singularity.nanocompore
	ln -s $(LAST_NANOCOMPORE_IMAGE) nanocompore_latest.img
