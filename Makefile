# Hermetic, verify-only build of the SQIsign library crates.
#
# This is a verification gate, not an artifact build: the multi-stage
# Containerfile builds and tests the crates offline in a digest-pinned stagex
# toolchain (`--locked`, `--network=none`), so a successful `docker build`
# means the workspace builds and its tests pass hermetically. There is no
# shipped image, so no OCI exporter is needed.
.PHONY: verify
verify:
	DOCKER_BUILDKIT=1 \
	docker build \
		--build-arg VERSION=$(VERSION) \
		--tag anchorageoss-sqisign-rs/sqisign-rs \
		--progress=plain \
		--label "org.opencontainers.image.source=https://github.com/anchorageoss/sqisign-rs" \
		$(if $(filter 1,$(NOCACHE)),--no-cache) \
		-f images/hermetic-build/Containerfile \
		.
