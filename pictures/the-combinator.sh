#!/bin/bash

# https://apple.stackexchange.com/questions/52879/how-to-combine-two-images-into-one-on-a-mac

# Append horizontally
convert +append figure-two-left-panel.jpg figure-two-right-panel.jpg 2016gl070016-p02.jpg

# Append vertically
# convert -append a.jpg b.jpg c.jpg
