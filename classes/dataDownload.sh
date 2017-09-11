wget -qO- -O tmp5.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FconvertIDs.db.zip?alt=media&token=55c80e8c-d5c1-43a5-995f-5e0c56242013' \
  && unzip tmp5.zip -d /srv/shiny-server/idep/data && rm tmp5.zip
wget -qO- -O tmp.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FgeneInfo.zip?alt=media&token=a281e5cf-6900-493c-81e4-89ce423c26bb'\
  && unzip tmp.zip -d /srv/shiny-server/idep/data && rm tmp.zip
wget -qO- -O tmp2.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fgmt.zip?alt=media&token=5ee100d1-e645-41ef-a591-7a9ba208ce3c' \
  && unzip tmp2.zip -d /srv/shiny-server/idep/data && rm tmp2.zip
wget -qO- -O tmp3.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2Fmotif.zip?alt=media&token=dc1e5972-ffd9-43a1-bbcc-49b78da6f047' \
  && unzip tmp3.zip -d /srv/shiny-server/idep/data && rm tmp3.zip
wget -qO- -O tmp4.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FpathwayDB.zip?alt=media&token=e602f2f7-102a-4cc4-8412-be2b05997daa' \
  && unzip tmp4.zip -d /srv/shiny-server/idep/data && rm tmp4.zip
wget -qO- -O tmp5.zip 'https://firebasestorage.googleapis.com/v0/b/firebase-bcloud.appspot.com/o/idep%2FgeneInfo%2FconvertIDs.db.zip?alt=media&token=55c80e8c-d5c1-43a5-995f-5e0c56242013' \
  && unzip tmp5.zip -d /srv/shiny-server/idep/data && rm tmp5.zip