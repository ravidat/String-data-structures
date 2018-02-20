 # command				output file
fmcount -i input-2.fm query-2 1 	output-2-1-error
fmcount -i input-2.fm query-2 2 	output-2
fmlocate -i input-2.fm query-2-locate 1	output-2-locate-1-error
fmlocate -i input-2.fm query-2-locate 2	output-2-locate


The time for "fmcount -i input-2.fm query-2 2" should be less than 4 seconds on the cicero1 server.
The size of input-2.fm should be ~800k.
