for file in j3.cpp j4.cpp p2.cpp p3.cpp p4.cpp s1.cpp s2.cpp s3.cpp s4.cpp t2.cpp t3.cpp t4.cpp ti2.cpp ti3.cpp ti4.cpp tr1.cpp tr2.cpp; do
    sed -i '' '1s/.*/#include "..\/..\/..\/..\/src\/join_ranked_results.cpp"/' "$file"
    sed -i '' '2s/.*//' "$file"
    sed -i '' '3s/.*//' "$file"
    sed -i '' '14s/.*//' "$file"
done

