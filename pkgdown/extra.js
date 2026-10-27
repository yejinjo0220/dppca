/*
 * Shorten only the function-name column of the pkgdown Reference index.
 * Do not change the underlying R help pages or function exports.
 * Works with pkgdown reference indexes rendered as <dl>/<dt> or <table>/<td>.
 */
(function () {
  "use strict";

  var displayNames = {
    "dp_pca": "dp_pca()",
    "plot.dp_pca": "plot.dp_pca()",
    "dp_loading": "dp_loading()",
    "plot.dp_loading": "plot.dp_loading()",
    "dp_scree": "dp_scree()",
    "plot.dp_scree": "plot.dp_scree()",
    "dp_score": "dp_score()",
    "plot.dp_score": "plot.dp_score()"
  };

  function getTopic(link) {
    try {
      var url = new URL(link.getAttribute("href"), document.baseURI);
      var match = url.pathname.match(/\/reference\/([^/]+)\.html$/);
      return match ? decodeURIComponent(match[1]) : null;
    } catch (e) {
      return null;
    }
  }

  function processEntry(entry) {
    var links = Array.from(entry.querySelectorAll("a[href]"));
    var chosen = links.find(function (link) {
      var topic = getTopic(link);
      return topic && Object.prototype.hasOwnProperty.call(displayNames, topic);
    });
    if (!chosen) return false;

    var topic = getTopic(chosen);
    // Keep the original click destination and its code styling.
    var code = chosen.closest("code");
    var element = code && entry.contains(code) ? code.cloneNode(true) : chosen.cloneNode(true);
    var anchor = element.tagName.toLowerCase() === "a" ? element : element.querySelector("a[href]");
    if (!anchor) return false;
    anchor.textContent = displayNames[topic];
    entry.replaceChildren(element);
    return true;
  }

  function cleanIndex() {
    // Limit this to the Reference landing page, not individual help pages.
    if (!/\/reference\/(?:index\.html)?$/.test(window.location.pathname)) return;
    var main = document.querySelector("main");
    if (!main) return;

    var entries = main.querySelectorAll(
      "dt, table tr > td:first-child, .ref-index .col-md-4, .ref-index .col-sm-4"
    );
    var changed = 0;
    entries.forEach(function (entry) {
      if (processEntry(entry)) changed += 1;
    });

    // A visible-in-DevTools diagnostic: should read "done:N" with N > 0.
    document.documentElement.setAttribute("data-dppca-ref-clean", "done:" + changed);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", cleanIndex);
  } else {
    cleanIndex();
  }
})();
