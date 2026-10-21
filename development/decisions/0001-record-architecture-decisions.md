# 0001. Record Architecture Decisions

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:** the user and Claude, following the convention shared with the laboratory and the sibling libraries

## 1. Context

physis is long-lived, and the laboratory depends on it. Future sessions need to know why its interfaces are as they are, including sessions that do not remember the original discussion.

## 2. Decision

- Record each decision that shapes the structure, interfaces or tooling of physis as a numbered record in `development/decisions/`.
- Never rewrite a record. A changed decision gets a new record that supersedes the old one.

## 3. Alternatives Considered

- **Commit messages only.** Rejected, because the reasons become hard to find.
- **A single design document edited in place.** Rejected, because it loses the history of why things changed.

## 4. Consequences

- Each decision costs a small amount of writing.
- `docs/architecture.md` describes the current state. The records explain how it came about.
