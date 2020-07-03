
# Build Docker image and push it to ECR

```
aws ecr get-login-password --region <region> | docker login --username AWS --password-stdin <ecr_repository>
docker build -t nanocompore/nanocompore_pipeline .
docker tag nanocompore/nanocompore_pipeline:latest <ecr_repository>/nanocompore/nanocompore_pipeline:latest
docker push <ecr_repository>/nanocompore/nanocompore_pipeline:latest
```

# Convert Docker image to Singularity

```
singularity build nanocompore_pipeline.img docker://nanocompore_pipeline:latest
```

