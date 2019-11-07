build() {
pandoc -t revealjs \
    --include-in-header simple-small.css \
    -V revealjs-url=https://revealjs.com \
    -V theme=white \
    -V transition=none \
    -o slides.html \
    --standalone slides.md
}

build

while :; do inotifywait -e modify slides.md; build; done

