resource "aws_launch_template" "nanocompore_launch_template" {
  name = "${var.pipeline_name}_lt"
  tag_specifications {
    resource_type = "instance"
    tags = {
      Project = var.pipeline_name
      ResourceGroup = var.resources_tag
    }
  }

 placement {
    availability_zone = "${var.region}${var.az}"
  }

  block_device_mappings {
    device_name = "/dev/xvda"

    ebs {
      volume_size = var.instance_volume_size
      volume_type = "gp2"
    }
  }

#  user_data = filebase64("${path.module}/efs_setup.sh")
   user_data = base64encode(<<EOF
MIME-Version: 1.0
Content-Type: multipart/mixed; boundary="==MYBOUNDARY=="

--==MYBOUNDARY==
Content-Type: text/x-shellscript; charset="us-ascii"

#!/bin/bash
yum install -y amazon-efs-utils wget unzip

wget "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip"
unzip awscli-exe-linux-x86_64.zip -d /
/aws/install --update -b /usr/bin

--==MYBOUNDARY==--

EOF
)
}

resource "aws_launch_template" "nanocompore_launch_template_fpga_nanopolish" {
  name = "${var.pipeline_name}_lt_fpga"
  tag_specifications {
    resource_type = "instance"
    tags = {
      Project = var.pipeline_name
      ResourceGroup = var.resources_tag
    }
  }

  block_device_mappings {
    device_name = "/dev/sda1"

    ebs {
      volume_size = var.instance_volume_size
      volume_type = "gp2"
    }
  }

#  user_data = filebase64("${path.module}/efs_setup.sh")
   user_data = base64encode(<<EOF
MIME-Version: 1.0
Content-Type: multipart/mixed; boundary="==MYBOUNDARY=="

--==MYBOUNDARY==
Content-Type: text/x-shellscript; charset="us-ascii"

#!/bin/bash
yum install -y amazon-efs-utils wget unzip

wget "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip"
unzip awscli-exe-linux-x86_64.zip -d /
/aws/install --update -b /usr/bin

--==MYBOUNDARY==--

EOF
)
}


