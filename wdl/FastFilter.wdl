version 1.0

task FastFilterSingleVcf {
  input {
    File input_vcf
    String output_vcf_name = "filtered.vcf.gz"
    String docker_image = "ghcr.io/broadinstitute/fastfilter:latest"
  }

  Int disk_size_gb = 1 + 3*ceil(size([input_vcf], "GB"))


  command <<<
    set -euo pipefail

    fastfilter filter \
      --input "~{input_vcf}" \
      --output "~{output_vcf_name}" \
      thresholds \
      --min-gq 5
  >>>

  output {
    File output_vcf = output_vcf_name
  }

  runtime {
    docker: docker_image
    memory: "8 GB"
    cpu: 4
    disks: "local-disk ~{disk_size_gb} SSD"
  }
}

workflow FastFilter {
  input {
    File input_vcf
    String output_vcf_name = "filtered.vcf.gz"
    String docker_image = "ghcr.io/broadinstitute/fastfilter:latest"
  }

  call FastFilterSingleVcf {
    input:
      input_vcf = input_vcf,
      output_vcf_name = output_vcf_name,
      docker_image = docker_image
  }

  output {
    File output_vcf = FastFilterSingleVcf.output_vcf
  }
}
