#!/bin/bash

# Script to backup/panTExp code to GitHub
# Usage: ./backup_to_github.sh [commit_message]
# Example: ./backup_to_github.sh "Added new feature"

REPO_PATH=$(pwd)
REMOTE_NAME="origin"
BRANCH_NAME="main"

# Function to get current timestamp
get_timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

# Get commit message from argument or use default
COMMENT="$1"
if [ -z "$COMMENT" ]; then
    COMMENT="Auto backup"
fi

echo "Starting backup to GitHub..."
echo "Repository: $REPO_PATH"
echo "Timestamp: $(get_timestamp)"
echo "Comment: $COMMENT"

# Navigate to repository directory
cd "$REPO_PATH" || exit 1

# Check if we're in a git repository
if [ ! -d ".git" ]; then
    echo "Error: Not a git repository"
    exit 1
fi

# Check for uncommitted changes
if [ -z "$(git status --porcelain)" ]; then
    echo "No changes to commit. Repository is up to date."
    exit 0
fi

# Add all changes
echo "Adding all changes..."
git add .

# Create commit with comment and timestamp
COMMIT_MSG="$COMMENT - $(get_timestamp)"
echo "Creating commit: $COMMIT_MSG"
git commit -m "$COMMIT_MSG"

# Push to remote
echo "Pushing to $REMOTE_NAME $BRANCH_NAME..."
git push "$REMOTE_NAME" "$BRANCH_NAME"

if [ $? -eq 0 ]; then
    echo "Backup completed successfully!"
else
    echo "Error: Backup failed"
    exit 1
fi
