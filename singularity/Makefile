LAST_COMMIT=$(shell git log -n 1 --pretty=format:%h -- Singularity.baseimage Singularity.nanocompore)
LAST_BASE_IMAGE=baseimage_$(LAST_COMMIT).img
LAST_NANOCOMPORE_IMAGE=nanocompore_$(LAST_COMMIT).img

all: nanocompore

.PHONY: nanocompore
nanocompore: $(LAST_NANOCOMPORE_IMAGE)
	
.PHONY: base
base: $(LAST_BASE_IMAGE)

$(LAST_BASE_IMAGE): Singularity.baseimage
	git diff-index --quiet HEAD -- Singularity.baseimage || ( echo "Commit the changes to Singularity.baseimage"; exit 1)
	sudo singularity build $(LAST_BASE_IMAGE) Singularity.baseimage
	ln -fs $(LAST_BASE_IMAGE) baseimage_latest.img

$(LAST_NANOCOMPORE_IMAGE): Singularity.nanocompore $(LAST_BASE_IMAGE)
	git diff-index --quiet HEAD -- Singularity.nanocompore || ( echo "Commit the changes to Singularity.nanocompore"; exit 1)
	sudo singularity build --writable $(LAST_NANOCOMPORE_IMAGE) Singularity.nanocompore
	ln -fs $(LAST_NANOCOMPORE_IMAGE) nanocompore_latest.img
