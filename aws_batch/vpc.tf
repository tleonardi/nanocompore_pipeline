resource "aws_vpc" "vpc" {
  cidr_block = "10.0.0.0/16"
  enable_dns_hostnames = "true"
}

resource "aws_subnet" "public_sn" {
  vpc_id     = aws_vpc.vpc.id
  map_public_ip_on_launch = true
  cidr_block = "10.0.1.0/24"
}

# Internet gateway for the public subnet
resource "aws_internet_gateway" "ig" {
  vpc_id = aws_vpc.vpc.id
  tags = {
    Name = "Nanocompore_pipeline_IG"
  }
}

resource "aws_security_group" "pipeline_sg" {
  name = "aws_batch_compute_environment_security_group"
  vpc_id = aws_vpc.vpc.id
  egress {
    from_port   = 0
    to_port     = 65535
    protocol    = "TCP"
    cidr_blocks = ["0.0.0.0/0"]
  }
  ingress {
    from_port   = 0
    to_port     = 65535
    protocol    = "TCP"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

# Routing table for public subnet
resource "aws_route_table" "public_sn_rt" {
  vpc_id = aws_vpc.vpc.id
  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = aws_internet_gateway.ig.id
  }
  tags = {
    Name = "public_sn_rt"
  }
}

# Associate the routing table to public subnet
resource "aws_route_table_association" "public_sn_rt_assn" {
  subnet_id = aws_subnet.public_sn.id
  route_table_id = aws_route_table.public_sn_rt.id
}

