# AWS Batch environment

This folder contains the Terraform code to automatically setup an AWS Batch environment suitable for running the Nanocompore pipeline.

The environment consists of three queues:

- a __main_queue__, which uses a mix of instance types depending on the rources required by the job
- a __gpu_queue__, which uses p3.2xlarge instance (equipped with one NVIDIA Tesla V100 GPU) for basecalling
- an __fpga_queue__, which uses f1.2xlarge instances (equipped with one Xilinx Virtex UltraScale+ VU9P FPGA) for hugenomic_nanopolish

## Usage

Rename `terraform.tfvars.template` to `terraform.tfvars` and configure the variables.

The following is a description of the variables defined in `terraform.tfvars`:

- `aws_access_key_id`: AWS access key
- `aws_secret_access_key`: AWS secret access key
- `aws_account_id`: AWS account ID 
- `bucket`: Name of the S3 bucket used by the pipeline 
- `pipeline_name`: Name of the pipeline 
- `instance_volume_size`: Size in GB of the EBS volume to attach to each instance 
- `resources_tag`: Value for the 'ResourceGroup' tag applied to all reasources created, including instances. 
- `region`: AWS region 
- `ami_id`: ID of the AMI to use 
- `gpu_ami_id`: ID of the GPU-optimised AMI to use 
- `fpga_ami_id`: ID of the FPGA-optimised AMI to use 
- `key_pair`: Name of the key pair to deploy. A keu under this name must be already available on AWS in the selected region 
 
## Choice of AMIs
The most recent ECS optimised AMIs can be found here: https://docs.aws.amazon.com/AmazonECS/latest/developerguide/ecs-optimized_AMI.html

Keep in mind that AMIs are region specific.

## Nanopolish with FPGA accelleration

[Huxelerate](https://www.huxelerate.it/) has developed an FPGA-accellerated version of Nanopolish which is available via a [dedicated AMI](https://aws.amazon.com/marketplace/pp/prodview-57t4mjkndu3sa).
