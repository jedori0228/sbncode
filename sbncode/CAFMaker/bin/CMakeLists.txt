cet_make_exec( extractCAFMetadata
               SOURCE extractCAFMetadata.cc
               LIBRARIES ${ROOT_BASIC_LIB_LIST}
               )

cet_script(diff_cafs)
cet_script(file_size_ana)

install_headers()
install_source()
