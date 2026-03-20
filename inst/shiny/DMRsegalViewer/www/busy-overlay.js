(function () {
  var overlayId = "dmrsegal-viewer-busy-overlay";
  var overlayActiveClass = "dmrsegal-viewer-overlay-active";
  var progressActiveClass = "dmrsegal-viewer-progress-active";
  var defaultMessage = "Processing...";
  var defaultDetail = "";
  var busyDelayMs = 250;

  function initBusyOverlay() {
    var overlay = document.getElementById(overlayId);
    if (!overlay || !document.body) {
      return;
    }

    var root = document.documentElement;
    var messageEl = overlay.querySelector('[data-role="message"]');
    var detailEl = overlay.querySelector('[data-role="detail"]');
    var cancelEl = overlay.querySelector('[data-role="cancel"]');
    var panelObserver = null;
    var bodyObserver = null;
    var busyTimer = null;
    var customState = { active: false };
    var customHandlerRegistered = false;

    function normalizeText(value) {
      return (value || "").replace(/\s+/g, " ").trim();
    }

    function textFrom(node, selector) {
      var target = node ? node.querySelector(selector) : null;
      return normalizeText(target ? target.textContent : "");
    }

    function isVisible(node) {
      if (!node || !node.isConnected) {
        return false;
      }

      var style = window.getComputedStyle(node);
      if (
        style.display === "none" ||
        style.visibility === "hidden" ||
        style.opacity === "0"
      ) {
        return false;
      }

      return node.getClientRects().length > 0;
    }

    function getProgressNotification() {
      var panel = document.getElementById("shiny-notification-panel");
      if (!panel) {
        return null;
      }

      var notifications = panel.querySelectorAll(".shiny-progress-notification");
      if (!notifications.length) {
        return null;
      }

      for (var i = notifications.length - 1; i >= 0; i -= 1) {
        if (isVisible(notifications[i])) {
          return notifications[i];
        }
      }

      return null;
    }

    function setOverlayText(message, detail) {
      if (messageEl) {
        messageEl.textContent = message || defaultMessage;
      }

      if (!detailEl) {
        return;
      }

      if (detail && detail.length) {
        detailEl.hidden = false;
        detailEl.textContent = detail;
      } else {
        detailEl.hidden = true;
        detailEl.textContent = "";
      }
    }

    function setCancelState(cancelable) {
      if (!cancelEl) {
        return;
      }

      cancelEl.hidden = !cancelable;
      cancelEl.disabled = !cancelable;
    }

    function showOverlay(message, detail, hasProgress, cancelable) {
      window.clearTimeout(busyTimer);
      busyTimer = null;

      setOverlayText(message, detail);
      setCancelState(!!cancelable);
      overlay.hidden = false;
      overlay.setAttribute("aria-hidden", "false");
      root.classList.add(overlayActiveClass);

      if (hasProgress) {
        root.classList.add(progressActiveClass);
      } else {
        root.classList.remove(progressActiveClass);
      }

      if (
        document.activeElement &&
        document.activeElement !== document.body &&
        typeof document.activeElement.blur === "function"
      ) {
        document.activeElement.blur();
      }
    }

    function hideOverlay() {
      window.clearTimeout(busyTimer);
      busyTimer = null;

      setCancelState(false);
      overlay.hidden = true;
      overlay.setAttribute("aria-hidden", "true");
      root.classList.remove(overlayActiveClass);
      root.classList.remove(progressActiveClass);
    }

    function applyCustomState() {
      if (!customState || !customState.active) {
        return false;
      }

      showOverlay(
        customState.message || defaultMessage,
        customState.detail || defaultDetail,
        false,
        !!customState.cancelable
      );
      return true;
    }

    function syncProgressOverlay() {
      if (applyCustomState()) {
        return true;
      }

      if (!root.classList.contains("shiny-busy")) {
        root.classList.remove(progressActiveClass);
        return false;
      }

      var progress = getProgressNotification();
      if (!progress) {
        root.classList.remove(progressActiveClass);
        return false;
      }

      showOverlay(
        textFrom(progress, ".progress-message") || defaultMessage,
        textFrom(progress, ".progress-detail"),
        true,
        false
      );
      return true;
    }

    function scheduleBusyFallback() {
      window.clearTimeout(busyTimer);
      busyTimer = window.setTimeout(function () {
        if (applyCustomState()) {
          return;
        }

        if (root.classList.contains("shiny-busy") && !getProgressNotification()) {
          showOverlay(defaultMessage, defaultDetail, false, false);
        }
      }, busyDelayMs);
    }

    function handleBusy() {
      if (applyCustomState()) {
        return;
      }

      if (!syncProgressOverlay()) {
        scheduleBusyFallback();
      }
    }

    function handleIdle() {
      if (applyCustomState()) {
        return;
      }

      if (!syncProgressOverlay()) {
        hideOverlay();
      }
    }

    function reconcileAfterOutput() {
      window.setTimeout(function () {
        if (applyCustomState()) {
          return;
        }

        if (root.classList.contains("shiny-busy")) {
          handleBusy();
        } else if (!syncProgressOverlay()) {
          hideOverlay();
        }
      }, 0);
    }

    function bindPanelObserver() {
      var panel = document.getElementById("shiny-notification-panel");

      if (panelObserver) {
        panelObserver.disconnect();
        panelObserver = null;
      }

      if (!panel) {
        if (!root.classList.contains("shiny-busy") && !applyCustomState()) {
          hideOverlay();
        }
        return;
      }

      panelObserver = new MutationObserver(function () {
        if (applyCustomState()) {
          return;
        }

        if (!syncProgressOverlay() && !root.classList.contains("shiny-busy")) {
          hideOverlay();
        }
      });

      panelObserver.observe(panel, {
        childList: true,
        subtree: true,
        characterData: true,
        attributes: true,
        attributeFilter: ["class", "style", "hidden", "aria-hidden"]
      });

      if (root.classList.contains("shiny-busy")) {
        handleBusy();
      } else if (!syncProgressOverlay()) {
        hideOverlay();
      }
    }

    function registerCustomHandler() {
      if (customHandlerRegistered) {
        return;
      }

      if (!window.Shiny || typeof window.Shiny.addCustomMessageHandler !== "function") {
        window.setTimeout(registerCustomHandler, 50);
        return;
      }

      window.Shiny.addCustomMessageHandler("dmrsegal-viewer-busy", function (message) {
        customState = message || { active: false };

        if (applyCustomState()) {
          return;
        }

        if (root.classList.contains("shiny-busy")) {
          handleBusy();
        } else if (!syncProgressOverlay()) {
          hideOverlay();
        }
      });

      customHandlerRegistered = true;
    }

    document.addEventListener("shiny:busy", handleBusy);
    document.addEventListener("shiny:idle", handleIdle);
    document.addEventListener("shiny:value", reconcileAfterOutput);
    document.addEventListener("shiny:error", reconcileAfterOutput);

    bodyObserver = new MutationObserver(function () {
      bindPanelObserver();
    });

    bodyObserver.observe(document.body, { childList: true });

    registerCustomHandler();
    bindPanelObserver();

    if (applyCustomState()) {
      return;
    }

    if (root.classList.contains("shiny-busy")) {
      handleBusy();
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initBusyOverlay, { once: true });
  } else {
    initBusyOverlay();
  }
})();
