#!/usr/bin/perl

# The script transposes a tab-delimited file.
# The script is called by the following bash script:
# transpose.smk

use strict;
use warnings;

# Initialize row and column counters
my $row_count = 0;
my $col_count = 0;
my $max_rows = 0;
my $max_cols = 0;
my @data = ();

# Read input line by line
while (<>) {
    chomp;  # Remove newline character
    my @fields = split /\t/;  # Split line into fields based on tab

    # Process each field in the current line
    for my $field (@fields) {
        if ($field =~ /./) {  # Check if the field is not empty
            $data[$row_count][$col_count] = $field;  # Store the field in the data array
        }
        $col_count++;  # Increment column counter
    }

    # Update maximum column count if necessary
    if ($col_count > $max_cols) {
        $max_cols = $col_count;
    }
    $col_count = 0;  # Reset column counter for the next row
    $row_count++;  # Increment row counter
}

$max_rows = $row_count;  # Set maximum row count

# Print the transposed data
for ($col_count = 0; $col_count < $max_cols; $col_count++) {
    for ($row_count = 0; $row_count < $max_rows; $row_count++) {
        if (defined($data[$row_count][$col_count])) {
            print $data[$row_count][$col_count];  # Print the transposed value
        }
        # Print a tab or newline based on the position
        if ($row_count == $max_rows - 1) {
            print "\n";  # Newline at the end of the row
        } else {
            print "\t";  # Tab between values
        }
    }
}