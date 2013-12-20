#!/bin/sh

#
# Read the item word list from performance_items.txt and 
# make performance_item.h.
# The list always start from etc and end with last
#


items=`cat performance_items.txt`

echo '//'
echo '// This is generated automatically.'
echo '// Edit performance_items.txt, not this file.'
echo '//'
echo '#ifndef PERFORMANCE_ITEMS_H'
echo '#define PERFORMANCE_ITEMS_H'
echo -n 'enum performance_item {etc, '

for item in $items
do
  echo $item | awk '{printf("%s, ", $1)}'
done
echo "last};"

echo -n 'static const char* performance_name[]= {"etc", '
for item in $items
do
  echo $item | awk '{printf("\"%s\", ", $1)}'
done

echo '"last"};'
echo '#endif'

