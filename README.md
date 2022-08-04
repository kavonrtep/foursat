# Tool for visualization of tandem repeats dotplots

## Step 1 - detection of similarities using blast converting blast output to matrix
```
usage: create_dotplot_matrix.py [-h] -i INPUT [-W WORD_LENGTHS] [-P PERC_IDENTITY] [-M MATRIX_WIDTH] [-S SCALE_COEF] [-o OUTPUT]

Run blast on sequence against itself and blast hits to dotplot matrix of similarities

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input fasta sequence
  -W WORD_LENGTHS, --word_lengths WORD_LENGTHS
                        input word size for blast separated by comma
  -P PERC_IDENTITY, --perc_identity PERC_IDENTITY
                        percentage identity thresholds separated by comma
  -M MATRIX_WIDTH, --matrix_width MATRIX_WIDTH
                        required dimension of output matrix
  -S SCALE_COEF, --scale_coef SCALE_COEF
                        Matrix_width is length_of_seq/scale_coef
  -o OUTPUT, --output OUTPUT
                        output directory, will be created

```


## Step 2 - plotting

```
plot_dotplot_matrix.py $input_directory
```

`$input_directory` is directory create by `create_dotplot_matrix.py` script.