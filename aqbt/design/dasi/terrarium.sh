CONTAINERID="4b940c5d726e"
FILEPATH="/Users/justinvrana/Desktop/test.dirt"
BINDTO="/data/test.dirt"
PASS="Mountain5"
USERNAME="vrana"
URL="http://52.27.43.242/"

docker run -v $FILEPATH:$BINDTO $CONTAINERID terrarium login --username $USERNAME --url $URL --password $PASS design --load-model-path "./docker/gibson.pkl" --filepath $BINDTO