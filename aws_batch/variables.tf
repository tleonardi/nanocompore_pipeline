# main creds for AWS connection
variable "aws_access_key_id" {
  description = "AWS access key"
}

variable "aws_secret_access_key" {
  description = "AWS secret access key"
}

variable "aws_account_id" {
  description = "AWS account ID"
}

variable "bucket" {
  description = "Name of the S3 bucket used by the pipeline"
}

variable "pipeline_name" {
  description = "Name of the pipeline"
}

variable "instance_volume_size" {
  description = "Size in GB of the EBS volume to attach to each instance"
} 

variable "resources_tag" {
  description = "Value for the 'ResourceGroup' tag applied to all reasources created, including instances."
}

variable "region" {
  description = "AWS region"
}

variable "ami_id" {
  description = "ID of the AMI to use"
}

variable "gpu_ami_id" {
  description = "ID of the GPU-optimised AMI to use"
}

variable "fpga_ami_id" {
  description = "ID of the FPGA-optimised AMI to use"
}

variable "key_pair" {
  description = "Name of the key pair to deploy. A keu under this name must be already available on AWS in the selected region"
}
