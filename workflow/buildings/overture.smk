"""
Download Overture Maps

Source
------


"""

rule download_overture:
  output:
    release_dir=directory("{OUTPUT_DIR}/input/overture/{params.RELEASE}")
  params:
    RELEASE="2024-01-17-alpha.0"
  shell:
    """
    pushd {output.release_dir}
      pwd
      # aws s3 sync \
      #   --region us-west-2 \
      #   --no-sign-request \
      #   s3://overturemaps-us-west-2/release/{params.RELEASE}/ .
    popd
    """
