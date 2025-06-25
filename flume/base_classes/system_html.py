def _write_html_file(output_directory, svg_filepath):
    """
    DOCS:
    """
    html_opener = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
    <meta charset="UTF-8">
    <title>Interactive Graphviz SVG</title>
    <script src="https://cdn.jsdelivr.net/npm/svg-pan-zoom@3.6.1/dist/svg-pan-zoom.min.js"></script>
    <style>
        .graph-node:hover polygon {
        stroke: blue;
        stroke-width: 2;
        cursor: pointer;
        }
    </style>
    </head>
    <body>

    <h2>Graph Visualization</h2>
    """

    with open(svg_filepath, "r", encoding="utf-8") as infile:
        svg_text = infile.read()

    # svg_html = f"""
    # <!-- Placeholder where SVG will be injected -->
    # <object type="image/svg+xml" data="{svg_filepath}">
    #   Your browser does not support SVG
    # </object>
    # """

    end_html = """
    <script>
    // Apply to all nodes with class 'graph-node'
    document.querySelectorAll('.graph-node').forEach(node => {
        node.addEventListener('click', (event) => {
        const title = node.querySelector('title')?.textContent || "Unknown node";
        alert(`You clicked on ${title}`);
        });
    });
    </script>

    <script>
    svgPanZoom('#my-svg', {
        zoomEnabled: true,
        controlIconsEnabled: true,
        fit: true,
        center: true,
        minZoom: 0.5,
        maxZoom: 10
    });
    </script>

    </body>
    </html>
    """

    # Write the html file
    with open(f"{output_directory}/test.html", "w") as f:
        f.write(html_opener)
        f.write(svg_text)
        f.write(end_html)

    print("HTML File created")

    return


def _write_html_file_test(output_directory, svg_filepath):
    html_content = """
        <!DOCTYPE html>
    <html lang="en">
    <head>
    <meta charset="UTF-8">
    <style>
        .graph-node:hover ellipse {
        stroke: red;
        stroke-width: 3;
        cursor: pointer;
        }
    </style>
    </head>
    <body>

    <!-- Inline or loaded SVG -->
    <svg id="graph-svg" xmlns="http://www.w3.org/2000/svg" width="500" height="500">
    <g id="node-A" class="graph-node">
        <title>A</title>
        <ellipse cx="100" cy="100" rx="50" ry="30" fill="lightblue"/>
        <text x="100" y="105" text-anchor="middle" fill="black">Node A</text>
    </g>
    <g id="node-B" class="graph-node">
        <title>B</title>
        <ellipse cx="300" cy="100" rx="50" ry="30" fill="lightgreen"/>
        <text x="300" y="105" text-anchor="middle" fill="black">Node B</text>
    </g>
    </svg>

    <script>
    // Apply to all nodes with class 'graph-node'
    document.querySelectorAll('.graph-node').forEach(node => {
        node.addEventListener('click', (event) => {
        const title = node.querySelector('title')?.textContent || "Unknown node";
        alert(`You clicked on title`);
        });
    });
    </script>

    </body>
    </html>
    """

    with open(f"{output_directory}/testv2.html", "w") as f:
        f.write(html_content)

    return
