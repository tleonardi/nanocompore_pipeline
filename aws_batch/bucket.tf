resource "aws_s3_bucket" "bucket" {
  bucket = var.bucket
  acl    = "private"

  tags = {
    Project = var.pipeline_name
    ResourceGroup = var.resources_tag
  }
}

