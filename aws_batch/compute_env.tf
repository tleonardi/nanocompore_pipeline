resource "aws_batch_compute_environment" "nanocompore_compute_env" {
  compute_environment_name_prefix = "${var.pipeline_name}_compute_env"
  tags = {
    Project = var.pipeline_name
    ResourceGroup = var.resources_tag
  }
  compute_resources {
    image_id = var.ami_id
    ec2_key_pair = var.key_pair
    instance_role = aws_iam_instance_profile.ecs_instance_role.arn

    max_vcpus = 128
    min_vcpus = 0
    instance_type = ["optimal"]

    security_group_ids = [
      aws_security_group.pipeline_sg.id,
    ]

    subnets = [
      aws_subnet.public_sn.id,
    ]

    type = "EC2"
    launch_template {
        launch_template_id = aws_launch_template.nanocompore_launch_template.id
        version = aws_launch_template.nanocompore_launch_template.latest_version
    }
  }

  service_role = aws_iam_role.aws_batch_service_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.aws_batch_service_role]
  lifecycle {
    create_before_destroy = true
  }
}

resource "aws_batch_job_queue" "main_queue" {
  name                 = "${var.pipeline_name}_main_q"
  state                = "ENABLED"
  priority             = 10
  compute_environments = [aws_batch_compute_environment.nanocompore_compute_env.arn]
  depends_on   = [aws_batch_compute_environment.nanocompore_compute_env]
}

resource "aws_batch_compute_environment" "nanocompore_compute_gpu_env" {
  compute_environment_name_prefix = "${var.pipeline_name}_compute_gpu_env"
  tags = {
    Project = var.pipeline_name
    ResourceGroup = var.resources_tag
  }

  compute_resources {
    image_id = var.gpu_ami_id
    ec2_key_pair = var.key_pair
    instance_role = aws_iam_instance_profile.ecs_instance_role.arn

    max_vcpus = 128
    min_vcpus = 0
    instance_type = ["p3.2xlarge"]

    security_group_ids = [
      aws_security_group.pipeline_sg.id,
    ]

    subnets = [
      aws_subnet.public_sn.id,
    ]

    type = "EC2"
    launch_template {
        launch_template_id = aws_launch_template.nanocompore_launch_template.id
        version = aws_launch_template.nanocompore_launch_template.latest_version
    }
  }

  service_role = aws_iam_role.aws_batch_service_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.aws_batch_service_role]
  lifecycle {
    create_before_destroy = true
  }
}

resource "aws_batch_job_queue" "gpu_queue" {
  name                 = "${var.pipeline_name}_gpu_q"
  state                = "ENABLED"
  priority             = 10
  compute_environments = [aws_batch_compute_environment.nanocompore_compute_gpu_env.arn]
  depends_on   = [aws_batch_compute_environment.nanocompore_compute_gpu_env]
}



resource "aws_batch_compute_environment" "nanocompore_compute_fpga_env" {
  compute_environment_name_prefix = "${var.pipeline_name}_compute_fpga_env"
  tags = {
    Project = var.pipeline_name
    ResourceGroup = var.resources_tag
  }

  compute_resources {
    image_id = var.fpga_ami_id
    ec2_key_pair = var.key_pair
    instance_role = aws_iam_instance_profile.ecs_instance_role.arn

    max_vcpus = 128
    min_vcpus = 0
    instance_type = ["f1.2xlarge"]

    security_group_ids = [
      aws_security_group.pipeline_sg.id,
    ]

    subnets = [
      aws_subnet.public_sn.id,
    ]

    type = "EC2"
    launch_template {
        launch_template_id = aws_launch_template.nanocompore_launch_template.id
        version = aws_launch_template.nanocompore_launch_template.latest_version
    }
  }

  service_role = aws_iam_role.aws_batch_service_role.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_role_policy_attachment.aws_batch_service_role]
  lifecycle {
    create_before_destroy = true
  }
}

resource "aws_batch_job_queue" "fpga_queue" {
  name                 = "${var.pipeline_name}_fpga_q"
  state                = "ENABLED"
  priority             = 10
  compute_environments = [aws_batch_compute_environment.nanocompore_compute_fpga_env.arn]
  depends_on   = [aws_batch_compute_environment.nanocompore_compute_fpga_env]
}

resource "aws_batch_job_definition" "nanopolish_hugenomic" {
  name = "nanopolish_hugenomic"
  type = "container"

  container_properties = <<CONTAINER_PROPERTIES
{
    "command": ["true"],
    "image": "tleonardi/nanopolish_hugenomic:latest",
    "memory": 4096,
    "vcpus": 8,
    "privileged": true,
    "user": "root",
    "mountPoints": [
       {
          "containerPath": "/usr/local/aws-cli/v2/current",
          "readOnly": true,
          "sourceVolume": "aws-cli"
       },
       {
          "containerPath": "/opt/xilinx",
          "readOnly": false,
          "sourceVolume": "vol-1"
       },
       {
          "containerPath": "/huxelerate",
          "readOnly": false,
          "sourceVolume": "vol-2"
       },
       {
          "containerPath": "/usr/share/hugenomic",
          "readOnly": false,
          "sourceVolume": "vol-3"
       },
       {
          "containerPath": "/usr/lib64",
          "readOnly": false,
          "sourceVolume": "vol-4"
       },
       {
          "containerPath": "/usr/lib",
          "readOnly": false,
          "sourceVolume": "vol-5"
       },
       {
          "containerPath": "/usr/bin",
          "readOnly": false,
          "sourceVolume": "vol-6"
       }
    ],
    "volumes": [
       {
          "host": {
             "sourcePath": "/usr/local/aws-cli/v2/current"
          },
          "name": "aws-cli"
       },
       {
          "host": {
             "sourcePath": "/opt/xilinx"
          },
          "name": "vol-1"
       },
       {
          "host": {
             "sourcePath": "/huxelerate"
          },
          "name": "vol-2"
       },
       {
          "host": {
             "sourcePath": "/usr/share/hugenomic"
          },
          "name": "vol-3"
       },
       {
          "host": {
             "sourcePath": "/usr/lib64"
          },
          "name": "vol-4"
       },
       {
          "host": {
             "sourcePath": "/usr/lib"
          },
          "name": "vol-5"
       },
       {
          "host": {
             "sourcePath": "/usr/bin"
          },
          "name": "vol-6"
       }
    ]
   }
CONTAINER_PROPERTIES
}


