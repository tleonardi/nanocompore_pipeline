provider "aws" {
  access_key  = var.aws_access_key_id
  secret_key  = var.aws_secret_access_key
  region      = var.region
}

resource "aws_iam_role" "ecs_instance_role" {
  name = "ecs_instance_role_${var.pipeline_name}"
  assume_role_policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
    {
        "Action": "sts:AssumeRole",
        "Effect": "Allow",
        "Principal": {
        "Service": "ec2.amazonaws.com"
        }
    }
    ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "ecs_instance_role_pol1" {
  role       = aws_iam_role.ecs_instance_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource "aws_iam_role_policy_attachment" "ecs_instance_role_pol2" {
  role       = aws_iam_role.ecs_instance_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonS3FullAccess"
}

resource "aws_iam_instance_profile" "ecs_instance_role" {
  name = "ecs_instance_role_${var.pipeline_name}"
  role = aws_iam_role.ecs_instance_role.name
}

resource "aws_iam_role" "aws_batch_service_role" {
  name = "aws_batch_service_role_${var.pipeline_name}"
  assume_role_policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
    {
        "Action": "sts:AssumeRole",
        "Effect": "Allow",
        "Principal": {
        "Service": "batch.amazonaws.com"
        }
    }
    ]
}
EOF
}

resource "aws_iam_role_policy_attachment" "aws_batch_service_role" {
  role       = aws_iam_role.aws_batch_service_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
}


resource "aws_iam_policy" "main_q_submission" {
  name        = "${var.pipeline_name}_main_q_exec_policy"
  path        = "/${var.pipeline_name}/"
  description = "Policy that grants job submission to the queues of pipeline ${var.pipeline_name}"

  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "0",
            "Effect": "Allow",
            "Action": [
                "batch:DeregisterJobDefinition",
                "batch:SubmitJob",
                "batch:RegisterJobDefinition"
            ],
            "Resource": [
                "arn:aws:batch:${var.region}:513365523772:job-queue/nanocompore_main_q",
                "arn:aws:batch:${var.region}:513365523772:job-queue/nanocompore_fpga_q",
                "arn:aws:batch:${var.region}:513365523772:job-queue/nanocompore_gpu_q",
                "arn:aws:batch:*:*:job-definition/*"
            ]
        },
        {
            "Sid": "1",
            "Effect": "Allow",
            "Action": [
                "batch:DescribeJobQueues",
                "batch:TerminateJob",
                "batch:DescribeJobs",
                "batch:CancelJob",
                "batch:DescribeJobDefinitions",
                "batch:ListJobs",
                "batch:DescribeComputeEnvironments"
            ],
            "Resource": "*"
        }
    ]
}
EOF
}


resource "aws_iam_policy" "bucket_access" {
  name        = "${var.pipeline_name}_bucket_access_policy"
  path        = "/${var.pipeline_name}/"
  description = "Policy that grants read/write access to the pipeline bucket"

  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:*"
       ],
      "Resource": [
        "${aws_s3_bucket.bucket.arn}",
        "${aws_s3_bucket.bucket.arn}/*"
      ]
    }
  ]
}
EOF
}

resource "aws_iam_group" "pipeline_users" {
  name = "${var.pipeline_name}_users"
  path = "/${var.pipeline_name}/"
}

resource "aws_iam_group_policy_attachment" "policy_for_bucket_access" {
  group      = aws_iam_group.pipeline_users.name
  policy_arn = aws_iam_policy.bucket_access.arn
}

resource "aws_iam_group_policy_attachment" "policy_for_job_submission" {
  group      = aws_iam_group.pipeline_users.name
  policy_arn = aws_iam_policy.main_q_submission.arn
}

