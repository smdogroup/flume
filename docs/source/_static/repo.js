document.addEventListener("DOMContentLoaded", function () {
  // Find the breadcrumb container
  const breadcrumbs = document.querySelector(".wy-breadcrumbs");
  if (breadcrumbs) {
    const repoLink = document.createElement("a");
    repoLink.href = "https://github.com/smdogroup/flume";
    repoLink.className = "repo-button";
    repoLink.target = "_blank";
    repoLink.rel = "noopener";
    repoLink.textContent = "View Repository";
    repoLink.style.marginLeft = "10px";
    breadcrumbs.appendChild(repoLink);
  }
});
