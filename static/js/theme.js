(function () {
  const root = document.documentElement;
  const toggle = document.getElementById("theme-toggle");
  if (!toggle) return;

  const stored = localStorage.getItem("theme");
  if (stored === "dark") {
    root.classList.add("dark");
  }

  toggle.addEventListener("click", function () {
    const isDark = root.classList.toggle("dark");
    localStorage.setItem("theme", isDark ? "dark" : "light");
  });
})();
