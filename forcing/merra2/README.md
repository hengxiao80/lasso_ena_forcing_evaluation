MERRA2 data can be downloaded using simple wget commands

The command to be used looks like this:

`wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies -i <url.txt>`

For this to work, you need a .netrc file and a .urs_cookies file in your home dir. The .netrc file contains your login and password to the NASA machine. The .urs_cookies file is just a cookie file that you can re-create (as empty using `touch .urs_cookies`) every time you do a new job.

<url.txt> contains a list of the files you want to download. This can be generated using the script `create_file_list.py`.

For LASSO ENA forcing generation for SAM, we need just 3 types of files: 
`MERRA2_400.inst3_3d_asm_Nv.*`,
`MERRA2_400.tavg1_2d_flx_Nx.*`,
`MERRA2_400.tavg1_2d_slv_Nx.*`.