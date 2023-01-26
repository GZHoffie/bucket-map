#!/bin/bash

# This script is used to find the correct number of buckets needed
# for bucketMap to set the correct template argument.
# Adopted from https://stackoverflow.com/questions/23992646/sequence-length-of-fasta-file.

# Get the name of the input file
FASTA_FILE=$1

# Maximum length of a bucket
BUCKET_LENGTH=$2


# Use awk to print the length of each sequence in the reference genome
awk -v bucket_l=${BUCKET_LENGTH} '
    function ceil(value)
    {
        return (value == int(value)) ? value : int(value)+1
    }

    /^>/ { # header pattern detected
        if (seqlen){
            # print previous seqlen if exists 
            buckets += ceil(seqlen/bucket_l)
        }

        # initialize sequence
        seqlen = 0

        # skip further processing
        next
      }

    # accumulate sequence length
    {
        seqlen += length($0)
    }

    # remnant seqlen if exists
    END{
        if(seqlen){
            buckets += ceil(seqlen/bucket_l)
        }
        print buckets
    }
' ${FASTA_FILE}