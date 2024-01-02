source setup_env.sh

ls $DST_SOURCE_PATH/$DST_NAME-000*-0000.root | cut -d '-' -f 2 | sed "s/^0*//" 
